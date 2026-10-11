# A raw query view is operation-owned. Retain reusable location state only
# for the actual published cache; detached raw views leave all session refs alone.
function _mesh_query_locator_locked(mesh::Union{Mesh,MixedMesh})
    mesh===LAST_MESH[] && return _cached_mesh_locator_locked(mesh)
    return mesh isa MixedMesh ? _MixedMeshLocator(mesh) : SimplexLocator(mesh)
end

# Build a detached view from explicit raw identity rows. Query readers receive
# its classification directly; they do not publish a generated mesh or alter
# node/element high-water marks, stored owners or original parameter records.
function _mesh_query_pair_locked(caller::AbstractString)
    model=_model_locked()
    cached=LAST_MESH[]
    if cached!==nothing
        class=_cached_classification_locked(cached)
        records=Model._discrete_mesh_records_model(model)
        any(!isempty(record.element_tags) for (_,record) in records) || return cached,class
        public=_mesh_query_cache_tags(cached,class)
        uncovered=any(any(tag->!haskey(public.element_indices,UInt64(tag)),record.element_tags)
                      for (_,record) in records)
        uncovered || return cached,class
        return _mesh_query_composed_pair(cached,class,public,records,model,caller)
    end
    records=Model._discrete_mesh_records_model(model)
    any((!isempty(record.node_tags) || !isempty(record.element_tags))
        for (_,record) in records) || throw(ArgumentError(
            "$caller: no mesh; call API.mesh.generate first"))
    staged=_reclassify_model_copy(model)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],ELEMENT_TAG_MAX[],typemax(Int32),typemax(Int32),String(caller))
    _dim01_ingest!(assembly,staged,nothing,nothing;preserve_source_order=true)
    plan=_dim01_finish!(assembly,staged,model_entities(staged),
        Set{Tuple{Int,Int32}}(),falses(length(assembly.nodes)),false,false;
        node_tag_order=true,record_view=true)
    class=_classification_with_public_tags(plan.class)
    return plan.mesh,class
end

# A legacy cache can coexist with caller-owned records. Its public label table
# is explicit even before the caller renames labels; overlay rows are actual
# nodes and therefore belong in the table's node count.
function _mesh_query_cache_tags(cached,class)
    class!==nothing && class.public_tags!==nothing && return class.public_tags
    overlay=_high_order_overlay(cached)
    count=overlay===nothing ? nnodes(cached) : nnodes(overlay)
    return _cache_public_tags(cached;node_count=count)
end

# This detached view carries unknown legacy owners in the private classification
# instead of inventing CAD entities in MixedEntityData. All node/cell identities
# are supplied by the old cache table or actual records. Coordinates never
# establish an identity. Readers leave every published object and high-water
# reference untouched.
function _mesh_query_composed_pair(cached,class,public,records,model,caller)
    actual,actual_class=_dim01_actual_cache(cached,class)
    points=NTuple{3,Float64}[Tuple(actual.coords[:,row]) for row in axes(actual.coords,2)]
    node_tags=copy(public.node_tags)
    node_positions=copy(public.node_indices)
    unknown=(4,Int32(0))
    stored=actual isa MixedMesh ? actual.entity_data : nothing
    node_owners=actual_class===nothing ?
                (stored===nothing ? fill(unknown,length(points)) : copy(stored.node_entities)) :
                copy(actual_class.node_entities)
    parameters=actual_class===nothing ?
        (stored===nothing ? Union{Nothing,Vector{Float64}}[nothing for _ in points] :
         deepcopy(stored.node_parametric)) :
        deepcopy(actual_class.node_parametric)
    for (entity,record) in sort!(collect(records);by=first)
        owner=(entity[1],Int32(entity[2]))
        for column in eachindex(record.node_tags)
            uid=UInt64(record.node_tags[column])
            point=ntuple(axis->record.node_coords[axis,column],3)
            position=get(node_positions,uid,0)
            parameter=Model._model_mesh_node_parameters(record,column)
            if position==0
                length(points)<typemax(Int32) || throw(ArgumentError(
                    "$caller: composed query exceeds the Int32 node limit"))
                all(isfinite,point) || throw(ArgumentError("$caller: source node $uid is non-finite"))
                push!(points,point);push!(node_tags,uid);push!(node_owners,owner)
                push!(parameters,parameter===nothing ? nothing : copy(parameter))
                node_positions[uid]=length(points)
            else
                points[position]==point || throw(ArgumentError(
                    "$caller: source node $uid conflicts with the published coordinates"))
                previous=node_owners[position]
                previous==unknown || previous==owner || throw(ArgumentError(
                    "$caller: source node $uid conflicts with the published owner"))
                parameters[position]===nothing && parameter!==nothing &&
                    (parameters[position]=copy(parameter))
            end
        end
    end
    order=Int[]
    columns=Dict{Int,Vector{Vector{Int32}}}()
    physical=Dict{Int,Vector{Int32}}()
    owners=Dict{Int,Vector{Int32}}()
    labels=Dict{Int,Vector{UInt64}}()
    cells=Dict{UInt64,Tuple{Int,Int32,Vector{Int32}}}()
    for (index,(msh,dim,offset,connectivity,classified)) in enumerate(_cache_catalog(actual,actual_class))
        isempty(connectivity) && continue
        kind=Int(msh)
        haskey(columns,kind) || push!(order,kind)
        group=get!(()->Vector{Int32}[],columns,kind)
        weights=get!(()->Int32[],physical,kind)
        entities=get!(()->Int32[],owners,kind)
        tags=get!(()->UInt64[],labels,kind)
        native_physical=_cache_native_blocks(actual)[index][3]
        for column in axes(connectivity,2)
            nodes=copy(connectivity[:,column])
            owner=actual_class!==nothing ? classified[column] : stored!==nothing ?
                stored.block_entities[index][column] :
                actual isa MixedMesh && actual.elementary_entities!==nothing ?
                    actual.elementary_entities[index][column] : Int32(0)
            uid=public.element_tags[offset+column]
            push!(group,nodes);push!(weights,native_physical[column])
            push!(entities,owner);push!(tags,uid)
            cells[uid]=(kind,owner,nodes)
        end
    end
    for (entity,record) in sort!(collect(records);by=first)
        owner=Int32(entity[2])
        for column in eachindex(record.element_tags)
            uid=UInt64(record.element_tags[column]);kind=Int(record.element_types[column])
            nodes=Int32[]
            for tag in record.element_nodes[column]
                position=get(node_positions,UInt64(tag),0)
                position>0 || throw(ArgumentError("$caller: source element $uid references unknown node $tag"))
                push!(nodes,Int32(position))
            end
            previous=get(cells,uid,nothing)
            if previous!==nothing
                previous[1]==kind && previous[3]==nodes &&
                    (previous[2]==0 || previous[2]==owner) || throw(ArgumentError(
                        "$caller: source element $uid conflicts with its published identity"))
                continue
            end
            haskey(columns,kind) || push!(order,kind)
            push!(get!(()->Vector{Int32}[],columns,kind),nodes)
            push!(get!(()->Int32[],physical,kind),_dim01_physical(model,entity))
            push!(get!(()->Int32[],owners,kind),owner)
            push!(get!(()->UInt64[],labels,kind),uid)
            cells[uid]=(kind,owner,nodes)
        end
    end
    coordinates=Matrix{Float64}(undef,3,length(points))
    for row in eachindex(points),axis in 1:3
        coordinates[axis,row]=points[row][axis]
    end
    blocks=ElementBlock[ElementBlock(kind,hcat(columns[kind]...),physical[kind]) for kind in order]
    element_tags=vcat((labels[kind] for kind in order)...)
    source=cached isa MixedMesh ? cached : nothing
    sections,partition=_mesh_query_composed_metadata(source,public,element_tags,caller)
    mesh=source===nothing ? MixedMesh(coordinates,blocks) :
        MixedMesh(coordinates,blocks;physical_names=source.physical_names,
            periodic_links=source.periodic_links,ancillary_sections=source.ancillary_sections,
            data_sections=sections,partition_data=partition)
    entities=Set{Tuple{Int,Int32}}()
    for kind in order,owner in owners[kind]
        owner>0 && push!(entities,(Int(msh_dimension(kind)),owner))
    end
    for owner in node_owners
        owner[1]<=3 && owner[2]>0 && push!(entities,owner)
    end
    boundaries=class===nothing ? Dict{Tuple{Int,Int32},Vector{Int32}}() : copy(class.boundaries)
    if class!==nothing
        union!(entities,class.entities)
    elseif stored!==nothing
        for entity in values(stored.entities)
            key=(entity.dim,entity.tag)
            push!(entities,key);boundaries[key]=copy(entity.boundaries)
        end
        if actual.partition_data!==nothing
            for entity in values(actual.partition_data.entities)
                key=(entity.dim,entity.tag)
                push!(entities,key);boundaries[key]=copy(entity.boundaries)
            end
        end
    end
    staged_class=_mixed_classification(mesh,(4,Int32(0)),sort!(collect(entities)),
        node_owners,boundaries,owners;node_parametric=parameters)
    table=_cache_public_tags(mesh;node_tags,element_tags,
                             authority=entities,caller)
    return mesh,_classification_with_public_tags(staged_class,table)
end

# Repeated-type query blocks are grouped, so preserved view-data and ghost-cell
# rows need an explicit old-public-UID to detached-row map. Cached node rows
# remain the prefix of the detached coordinate array and retain their indices.
function _mesh_query_composed_metadata(source,public,element_tags,caller)
    source===nothing && return Elements.MshDataSection[],nothing
    (isempty(source.data_sections) && source.partition_data===nothing) &&
        return Elements.MshDataSection[],nothing
    positions=Dict(uid=>Int32(row) for (row,uid) in enumerate(element_tags))
    cell_map=Int32[positions[uid] for uid in public.element_tags]
    sections=isempty(source.data_sections) ? Elements.MshDataSection[] :
        _reclassify_remap_sections(source,Int32.(1:nnodes(source)),cell_map,caller)
    previous=source.partition_data
    partition=if previous===nothing
        nothing
    else
        ghosts=Elements.MixedGhostElement[Elements.MixedGhostElement(cell_map[ghost.element],
            ghost.partition,ghost.ghost_partitions) for ghost in previous.ghost_elements]
        Elements.MixedPartitionData(previous.num_partitions;entities=previous.entities,
            ghost_entities=previous.ghost_entities,ghost_elements=ghosts)
    end
    return sections,partition
end

_mesh_query_detached_public(mesh,class)=mesh!==LAST_MESH[] && mesh isa MixedMesh &&
    mesh.entity_data===nothing && class!==nothing ? class.public_tags : nothing

# Global edge/face identifiers refer to the published cache's row space. Rebind
# owned copies to the detached view by public node identity before extending a
# hierarchical-key catalogue; queries never replace the published catalogues.
function _mesh_query_key_catalogs(mesh,class)
    edges=LAST_MESH_EDGES[];faces=LAST_MESH_FACES[]
    mesh===LAST_MESH[] && return edges,faces
    (edges===nothing && faces===nothing) && return nothing,nothing
    source=LAST_MESH[]
    source===nothing && return nothing,nothing
    public=class===nothing ? nothing : class.public_tags
    public===nothing && throw(ArgumentError("API.mesh key query: detached view has no public node identities"))
    original=_mesh_query_cache_tags(source,_cached_classification_locked(source))
    rows=Int32[public.node_indices[tag] for tag in original.node_tags]
    count=nnodes(mesh)
    edge_copy=if edges===nothing
        nothing
    else
        nodes=NTuple{2,Int32}[ntuple(slot->rows[node[slot]],2) for node in edges.nodes]
        tags=Dict(MeshEntityTopology._edge_key(node...)=>identifier
            for (node,identifier) in zip(nodes,edges.identifiers))
        MeshEntityTopology.MeshEdgeTopology(count,nodes,copy(edges.identifiers),tags,copy(edges.used_tags))
    end
    face_copy=if faces===nothing
        nothing
    else
        triangles=NTuple{3,Int32}[ntuple(slot->rows[node[slot]],3) for node in faces.triangle_nodes]
        quadrangles=NTuple{4,Int32}[ntuple(slot->rows[node[slot]],4) for node in faces.quadrangle_nodes]
        triangle_tags=Dict(MeshEntityTopology._triangle_key(node...)=>identifier
            for (node,identifier) in zip(triangles,faces.triangle_identifiers))
        quadrangle_tags=Dict(MeshEntityTopology._quadrangle_key(node...)=>identifier
            for (node,identifier) in zip(quadrangles,faces.quadrangle_identifiers))
        MeshEntityTopology.MeshFaceTopology(count,triangles,copy(faces.triangle_identifiers),
            triangle_tags,quadrangles,copy(faces.quadrangle_identifiers),quadrangle_tags,
            copy(faces.used_tags))
    end
    return edge_copy,face_copy
end

# Adding independent records must not replace a legacy cell's numerical
# provider. The detached catalogue only supplies public identities and column
# order; retained simplex cells still use their original Mesh or P2 overlay.
function _mesh_query_legacy_block(source::Mesh,msh,caller)
    overlay=_high_order_overlay(source)
    if overlay!==nothing
        msh==_p2_etype(overlay) && return overlay,_p2_tag_offset(source,overlay),_p2_cells(overlay)
        msh==_p2_skeleton_etype(overlay) && return nothing
    end
    block=_mesh_element_block(source,msh)
    return block===nothing ? nothing : (source,block[1],block[2])
end

function _mesh_query_scatter_jacobians!(target,values,slots,point_count)
    for component in 1:3
        width=component==1 ? 9 : component==2 ? 1 : 3
        stride=width*point_count
        for (position,slot) in enumerate(slots),entry in 1:stride
            target[component][stride*(slot-1)+entry]=values[component][stride*(position-1)+entry]
        end
    end
    return nothing
end

function _mesh_query_composed_jacobians(mesh,class,msh,local_coord,positions,caller)
    source=LAST_MESH[]
    (source isa Mesh && mesh!==source && mesh isa MixedMesh) || return nothing
    table=class===nothing ? nothing : class.public_tags
    table===nothing && return nothing
    legacy=_mesh_query_legacy_block(source,msh,caller)
    legacy===nothing && return nothing
    block=_mesh_element_block(mesh,msh)
    block===nothing && return nothing
    offset,cells=block
    local_coordinates,point_count=MeshReferenceGeometry._checked_local_coordinates(local_coord,caller)
    selected=positions===nothing ? (1:size(cells,2)) : positions isa UnitRange{Int} ?
        (MeshReferenceGeometry._checked_element_range(positions,size(cells,2),caller);positions) :
        MeshReferenceGeometry._checked_element_positions(positions,size(cells,2),caller)
    (isempty(selected) || point_count==0) && return (Float64[],Float64[],Float64[])
    original=_mesh_query_cache_tags(source,_cached_classification_locked(source))
    geometry,legacy_offset,legacy_cells=legacy
    old_columns=Int[];old_slots=Int[];raw_columns=Int[];raw_slots=Int[]
    for (slot,column) in enumerate(selected)
        uid=table.element_tags[offset+column]
        dense=get(original.element_indices,uid,0)
        if legacy_offset<dense<=legacy_offset+size(legacy_cells,2)
            push!(old_columns,dense-legacy_offset);push!(old_slots,slot)
        else
            push!(raw_columns,column);push!(raw_slots,slot)
        end
    end
    isempty(old_columns) && return nothing
    result=MeshReferenceGeometry._allocate_results(length(selected),point_count,caller)
    old_values=geometry isa Mesh ? mesh_jacobians(geometry,msh,local_coordinates,old_columns) :
        _legacy_p2_jacobians(geometry,local_coordinates,old_columns,legacy_offset,caller)
    _mesh_query_scatter_jacobians!(result,old_values,old_slots,point_count)
    if !isempty(raw_columns)
        raw_values=mesh_jacobians(mesh,msh,local_coordinates,raw_columns)
        _mesh_query_scatter_jacobians!(result,raw_values,raw_slots,point_count)
    end
    return result
end

function _mesh_query_composed_jacobian(mesh,class,tag,local_coord,caller)
    (LAST_MESH[] isa Mesh && mesh!==LAST_MESH[] && mesh isa MixedMesh) || return nothing
    msh,_,_,cell,_=_mixed_query_reference(mesh,tag,caller)
    return _mesh_query_composed_jacobians(mesh,class,msh,local_coord,cell:cell,caller)
end

# GEdge bulk discovery advertises the first actual Line factory in an entity.
# Other stored Line factories retain their identity and remain individually
# queryable; they are omitted only from this native bulk discovery operation.
function _mesh_curve_first_type(model,cached,class,entity;ghost_positions=nothing)
    if cached!==nothing && class!==nothing
        receipt=get(model.meshing.imported_entities,(1,Int(entity)),nothing)
        ghost=receipt!==nothing && startswith(receipt.kind,"Ghost ")
        for (msh,dim,offset,cells,owners) in _cache_catalog(cached,class)
            dim==1 || continue
            if ghost
                isempty(_public_bulk_cell_positions(model,dim,owners,entity;offset,ghost_positions)) || return Int32(msh)
            elseif any(==(Int32(entity)),owners)
                return Int32(msh)
            end
        end
    end
    record=_dim01_record(model,(1,Int(entity)))
    return record===nothing || isempty(record.element_types) ? nothing :
        first(record.element_types)
end

function _record_bulk_element_blocks(record,dimension)
    types,tags,nodes=_record_element_blocks(record)
    dimension==1 && !isempty(record.element_types) || return types,tags,nodes
    position=findfirst(==(first(record.element_types)),types)
    position===nothing && return Int32[],Vector{UInt64}[],Vector{UInt64}[]
    return types[position:position],tags[position:position],nodes[position:position]
end
