# The emitted constructor receipt authenticates a detached geometry view.
# Publish a fresh receipt after the final authority and emitted parts exist.
# Each later authority replacement needs a new receipt; an old receipt cannot
# authenticate a rebuilt/reclassified/current generation by shape alone.
struct _ModelMeshEmittedReceipt
    actual::MixedMesh
    carrier::Union{Mesh,MixedMesh}
    parts::Dict{Tuple{Int,Int32},Tuple{Union{Mesh,MixedMesh},Vector{Int32}}}
    cell_owners::Dict{Tuple{Int,Tuple},Tuple{Int,Int32}}
    carrier_rows::Vector{Int32}
    native_state::NamedTuple
end

_emitted_projection_blocks(mesh::Mesh)=((1,mesh.segs,mesh.seg_tag),
    (2,mesh.tris,mesh.tri_tag),(4,mesh.tets,mesh.tet_tag))
_emitted_projection_blocks(mesh::MixedMesh)=
    ((Int(block.msh),block.nodes,block.tags) for block in mesh.blocks)

function _model_mesh_emitted_receipt(authority,parts,part_rows,caller)
    authority.mesh!==nothing && authority.carrier!==nothing || throw(ArgumentError(
        "$caller: emitted receipt requires an actual authority and its carrier"))
    length(parts)==length(part_rows) || throw(ArgumentError(
        "$caller: emitted constructor parts do not align with their row bindings"))
    bound=Dict{Tuple{Int,Int32},Tuple{Union{Mesh,MixedMesh},Vector{Int32}}}()
    owners=Dict{Tuple{Int,Tuple},Tuple{Int,Int32}}()
    for ((dimension,tag,part),rows) in zip(parts,part_rows)
        owner=(dimension,Int32(tag))
        haskey(authority.products,owner) && !haskey(bound,owner) || throw(ArgumentError(
            "$caller: emitted constructor part has absent or duplicate ownership"))
        length(rows)==nnodes(part) && allunique(rows) &&
            all(row->0<row<=nnodes(authority.carrier),rows) || throw(ArgumentError(
            "$caller: emitted constructor rows lack a unique carrier binding"))
        isequal(part.coords,authority.carrier.coords[:,rows]) || throw(ArgumentError(
            "$caller: emitted constructor coordinates disagree with their bound carrier"))
        part isa MixedMesh && part.entity_data!==nothing && throw(ArgumentError(
            "$caller: an emitted constructor part claims actual native metadata"))
        all(entry->isempty(entry[2]) || msh_dimension(entry[1])==dimension,
            _emitted_projection_blocks(part)) || throw(ArgumentError(
            "$caller: emitted constructor part contains a foreign cell dimension"))
        for (msh,cells,tags) in _emitted_projection_blocks(part),column in axes(cells,2)
            all(row->0<row<=length(rows),cells[:,column]) || throw(ArgumentError(
                "$caller: emitted constructor cell references an absent bound vertex"))
            key=(msh,Tuple(rows[cells[:,column]]))
            haskey(owners,key) && throw(ArgumentError("$caller: emitted constructors repeat a typed carrier cell"))
            owners[key]=owner
        end
        bound[owner]=(deepcopy(part),Int32.(rows))
    end
    carrier_keys=Set((msh,Tuple(cells[:,column])) for (msh,cells,tags) in
        _emitted_projection_blocks(authority.carrier) for column in axes(cells,2))
    Set(keys(owners))==carrier_keys || throw(ArgumentError(
        "$caller: emitted constructors do not partition the current carrier"))
    return _ModelMeshEmittedReceipt(authority.mesh,authority.carrier,bound,owners,
        copy(authority.carrier_rows),_model_mesh_emitted_state(authority.mesh))
end

# Snapshot only native state that authenticates row/cell identity. Coordinate
# arithmetic of an emitted view never rewrites this actual state.
function _model_mesh_emitted_state(mesh)
    data=mesh.entity_data
    return (coords=copy(mesh.coords),node_entities=copy(data.node_entities),
        node_parametric=deepcopy(data.node_parametric),node_tags=copy(data.external_node_tags),
        entities=deepcopy(data.entities),links=deepcopy(mesh.periodic_links),
        cells=[(Int(block.msh),copy(block.nodes),copy(block.tags),copy(owners),copy(tags))
            for (block,owners,tags) in zip(mesh.blocks,data.block_entities,data.external_element_tags)])
end

function _model_mesh_emitted_state_matches(mesh,state)
    data=mesh.entity_data
    isequal(mesh.coords,state.coords) && isequal(data.node_entities,state.node_entities) &&
        isequal(data.node_parametric,state.node_parametric) &&
        isequal(data.external_node_tags,state.node_tags) &&
        Set(keys(data.entities))==Set(keys(state.entities)) || return false
    for (owner,entity) in data.entities,field in fieldnames(MixedEntity)
        isequal(getfield(entity,field),getfield(state.entities[owner],field)) || return false
    end
    length(mesh.periodic_links)==length(state.links) || return false
    for (left,right) in zip(mesh.periodic_links,state.links),field in fieldnames(MixedPeriodicLink)
        isequal(getfield(left,field),getfield(right,field)) || return false
    end
    current=[(Int(block.msh),block.nodes,block.tags,owners,tags)
        for (block,owners,tags) in zip(mesh.blocks,data.block_entities,data.external_element_tags)]
    return isequal(current,state.cells)
end

function _model_mesh_with_emitted_receipt(authority,parts,rows;caller="emitted geometry")
    receipt=_model_mesh_emitted_receipt(authority,parts,rows,caller)
    return _ModelMeshAuthority(authority.mesh,authority.products,authority.done,
        authority.private_rank,authority.native_inputs,authority.carrier,authority.carrier_rows,
        authority.producer_rows,authority.carrier_cells,authority.surface_constructor_inputs,receipt)
end

# These branches still require the current generation's exact constructor
# receipt and literal cell incidence below. Other naked inputs retain their
# existing geometric certification path.
function _model_mesh_emitted_scope(m,tag)
    params=_extrude_gate(m,3,tag)
    if params!==nothing
        return params.recombine && params.quad_to_tri===:none
    end
    haskey(m.volumes,tag) && !isempty(m.volumes[tag]) &&
        !haskey(m.meshing.transfinite_volumes,tag) && !_implicit_volume_surface(m,tag)
end

function _model_mesh_emitted_projection(m,mesh,dimension,tag,receipt,caller)
    receipt isa _ModelMeshEmittedReceipt || return nothing
    dimension==3 || return nothing
    _model_mesh_emitted_scope(m,tag) || return nothing
    mesh isa MixedMesh && mesh.entity_data!==nothing && return nothing
    authority=m.meshing.mesh_authority
    authority!==nothing && authority.mesh===receipt.actual &&
        authority.carrier===receipt.carrier || return nothing
    owner=(dimension,Int32(tag))
    constructor=get(receipt.parts,owner,nothing)
    constructor===nothing && return nothing
    part,constructor_rows=constructor
    # Exact constructor equality authenticates the detached naked part.
    # Coordinates verify those row events; they do not supply native identity.
    nnodes(mesh)==nnodes(part) || return nothing
    supplied=collect(_emitted_projection_blocks(mesh))
    expected=collect(_emitted_projection_blocks(part))
    length(supplied)==length(expected) && all(isequal(a,b) for (a,b) in zip(supplied,expected)) ||
        return nothing
    isequal(mesh.coords,part.coords) || throw(ArgumentError(
        "$caller: emitted constructor coordinates no longer match their selected row receipt"))
    isequal(mesh.coords,authority.carrier.coords[:,constructor_rows]) || throw(ArgumentError(
        "$caller: emitted constructor receipt no longer matches its carrier coordinates"))
    authority.carrier_rows==receipt.carrier_rows &&
        _model_mesh_emitted_state_matches(authority.mesh,receipt.native_state) || throw(ArgumentError(
        "$caller: emitted receipt no longer identifies its unchanged actual constructor state"))
    for (owned_part,rows) in values(receipt.parts)
        isequal(owned_part.coords,authority.carrier.coords[:,rows]) || throw(ArgumentError(
            "$caller: emitted closure coordinates no longer match their constructor row receipt"))
    end
    actual=authority.mesh;data=actual.entity_data
    closure=_model_mesh_projection_closure(m,dimension,tag,caller)
    actual_rows=authority.carrier_rows[constructor_rows]
    allunique(actual_rows) && Set(actual_rows)==Set(authority.products[owner].rows) ||
        throw(ArgumentError("$caller: emitted part does not bind the complete selected actual entity"))
    # Verify actual cell UID, owner, family and ordered primary incidence for
    # every usable carrier cell. Proxy/subdivision cells need their own planner.
    native=Dict{UInt64,Tuple{Int,Tuple{Int,Int32},Vector{Int32}}}()
    for (block,owners,uids) in zip(actual.blocks,data.block_entities,data.external_element_tags)
        for column in axes(block.nodes,2)
            uid=uids[column]
            haskey(native,uid) && throw(ArgumentError("$caller: duplicate actual cell UID"))
            native[uid]=(Int(block.msh),(msh_dimension(block.msh),owners[column]),copy(block.nodes[:,column]))
        end
    end
    carrier_cells=collect(_emitted_projection_blocks(authority.carrier))
    length(carrier_cells)==length(authority.carrier_cells) || throw(ArgumentError(
        "$caller: emitted carrier cell blocks lack their actual constructor events"))
    events=Dict{Tuple{Int,Tuple},UInt64}();seen=Set{UInt64}()
    for ((msh,cells,tags),uids) in zip(carrier_cells,authority.carrier_cells)
        length(uids)==size(cells,2) || throw(ArgumentError(
            "$caller: emitted carrier cell columns lack their actual constructor events"))
        for column in axes(cells,2)
            key=(msh,Tuple(cells[:,column]))
            cell_owner=get(receipt.cell_owners,key,nothing)
            cell_owner!==nothing || throw(ArgumentError("$caller: carrier cell lacks its emitted constructor owner"))
            cell_owner in closure || continue
            uid=uids[column]
            uid!==nothing || throw(ArgumentError("$caller: emitted closure cell lacks its actual constructor UID"))
            current=get(native,uid,nothing)
            current!==nothing && current[1]==msh && current[2]==cell_owner &&
                current[3]==authority.carrier_rows[cells[:,column]] || throw(ArgumentError(
                "$caller: emitted carrier cell loses its typed actual UID incidence"))
            uid in seen && throw(ArgumentError("$caller: emitted carrier repeats an actual cell event"))
            push!(seen,uid)
            key=(msh,Tuple(cells[:,column]))
            haskey(events,key) && throw(ArgumentError("$caller: emitted carrier repeats a typed constructor cell"))
            events[key]=uid
        end
    end
    selected=Set{UInt64}()
    for (msh,cells,tags) in supplied,column in axes(cells,2)
        uid=get(events,(msh,Tuple(constructor_rows[cells[:,column]])),nothing)
        uid!==nothing && native[uid][2]==owner || throw(ArgumentError(
            "$caller: emitted selected cell lacks its unique owned actual constructor event"))
        push!(selected,uid)
    end
    selected==Set(uid for (uid,cell) in native if cell[2]==owner) || throw(ArgumentError(
        "$caller: emitted selected cells do not cover the current actual entity"))
    dimension==3 || throw(ArgumentError("$caller: this emitted receipt adapter requires the Volume planner"))
    _volume_projection_contract(mesh,caller)
    any(link->(link.dim,link.slave_entity) in closure || (link.dim,link.master_entity) in closure,
        actual.periodic_links) && throw(ArgumentError(
        "$caller: emitted periodic closure requires its affine constructor planner"))
    inverse=Dict(row=>Int32(position) for (position,row) in enumerate(authority.carrier_rows))
    rows=copy(actual_rows);present=Set(rows)
    needed=Set{Int32}()
    for selected_owner in closure
        product=get(authority.products,selected_owner,nothing)
        product===nothing || union!(needed,product.rows)
    end
    for row in sort!(collect(setdiff(needed,present));by=row->authority.private_rank[row])
        haskey(inverse,row) || throw(ArgumentError(
            "$caller: emitted closure loses an actual owned vertex row"))
        push!(rows,row)
    end
    row_map=Dict(row=>Int32(position) for (position,row) in enumerate(rows))
    coords=authority.carrier.coords[:,Int32[inverse[row] for row in rows]]
    blocks=ElementBlock[];owners=Vector{Int32}[]
    for (block,entities,uids) in zip(actual.blocks,data.block_entities,data.external_element_tags)
        columns=Int[column for column in axes(block.nodes,2) if
            msh_dimension(block.msh)<dimension && (msh_dimension(block.msh),entities[column]) in closure]
        isempty(columns) && continue
        for column in columns
            # Constructor-UID attestation is only possible for cells whose
            # owner entity was itself emitted as a part: emitted parts
            # partition the carrier, so their lower cells appear in `seen`.
            # Un-emitted lower entities (and Point cells, which are never
            # carrier cells) are certified by the bound-row, typed-boundary
            # and edge-complex checks below instead.
            cell_owner=(msh_dimension(block.msh),entities[column])
            (cell_owner[1]==0 || !haskey(receipt.parts,cell_owner) ||
                uids[column] in seen) || throw(ArgumentError(
                "$caller: emitted lower cell lacks its literal actual constructor event"))
        end
        nodes=map(row->get(row_map,row,Int32(0)),block.nodes[:,columns])
        all(>(0),nodes) || throw(ArgumentError("$caller: emitted lower cell references an absent bound row"))
        physical=Int32[_model_projection_legacy_tag(_model_projection_physical_tags(m,
            msh_dimension(block.msh),Int(entities[column]))) for column in columns]
        push!(blocks,ElementBlock(block.msh,nodes,physical));push!(owners,entities[columns])
    end
    physical=_model_projection_legacy_tag(_model_projection_physical_tags(m,dimension,tag))
    for (msh,cells,tags) in supplied
        isempty(cells) && continue
        push!(blocks,ElementBlock(msh,cells,fill(physical,size(cells,2))))
        push!(owners,fill(Int32(tag),size(cells,2)))
    end
    # The projection still requires the complete actual typed boundary and
    # every literal lower Line in the emitted Volume's edge complex.
    boundary_surfaces=Set(Int32.(abs.(_model_volume_boundary_surfaces(m,tag,caller))))
    lower_faces=Set{Tuple}();line_edges=Dict{Tuple,Int32}()
    for (block,block_owners) in zip(blocks,owners),column in axes(block.nodes,2)
        cell_dimension=msh_dimension(block.msh)
        if cell_dimension==2 && block_owners[column] in boundary_surfaces
            face=Tuple(sort!(collect(block.nodes[:,column])))
            face in lower_faces && throw(ArgumentError("$caller: emitted lower Surfaces repeat a typed boundary facet"))
            push!(lower_faces,face)
        elseif cell_dimension==1
            edge=Tuple(sort!(collect(block.nodes[:,column])))
            previous=get(line_edges,edge,nothing)
            previous===nothing || previous==block_owners[column] || throw(ArgumentError(
                "$caller: emitted lower Line incidence aliases distinct Curve owners"))
            line_edges[edge]=block_owners[column]
        end
    end
    lower_faces==Set(Tuple(face) for face in first(_volume_cell_boundary_faces(mesh))) ||
        throw(ArgumentError("$caller: emitted lower Surfaces do not cover its exact typed Volume boundary"))
    volume_edges=_volume_edge_set(mesh)
    all(edge->edge in volume_edges,keys(line_edges)) || throw(ArgumentError(
        "$caller: emitted lower Line is absent from its Volume edge complex"))
    # Synthetic projection labels preserve the geometric frontend contract;
    # actual UIDs served solely to attest the constructor and lower incidence.
    entities=Dict{Tuple{Int,Int},MixedEntity}();names=Dict{Tuple{Int,Int},String}()
    for selected_owner in closure
        current=get(data.entities,selected_owner,nothing)
        current!==nothing || throw(ArgumentError("$caller: emitted closure lacks entity metadata"))
        entity_nodes=Set{Int32}(Int32(position) for (position,row) in enumerate(rows)
            if data.node_entities[row]==selected_owner)
        for (block,block_owners) in zip(blocks,owners),column in axes(block.nodes,2)
            (msh_dimension(block.msh),block_owners[column])==selected_owner || continue
            union!(entity_nodes,block.nodes[:,column])
        end
        tags=_model_projection_physical_tags(m,selected_owner[1],Int(selected_owner[2]))
        entities[selected_owner]=_model_mesh_entity_metadata(m,selected_owner[1],Int(selected_owner[2]),
            _model_projection_bbox((;coords),entity_nodes,caller);physical_tags=tags,
            boundaries=copy(current.boundaries),embedded_curves=copy(current.embedded_curves))
        for physical_tag in tags
            key=(selected_owner[1],Int(physical_tag))
            haskey(m.physical_names,key) && (names[key]=m.physical_names[key])
        end
    end
    metadata=MixedEntityData(entities;
        node_entities=data.node_entities[rows],node_parametric=Union{Nothing,Vector{Float64}}[nothing for _ in rows],
        external_node_tags=UInt64.(1:length(rows)),block_entities=owners,
        external_element_tags=_model_projection_external_elements(blocks))
    output=MixedMesh(coords,blocks;physical_names=names,entity_data=metadata,elementary_entities=owners)
    diagnostic=validate(output)
    diagnostic.ok || throw(ArgumentError("$caller: invalid emitted projection — "*join(diagnostic.messages,"; ")))
    return output
end
