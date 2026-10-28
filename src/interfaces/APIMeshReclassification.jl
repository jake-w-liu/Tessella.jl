# Records are read throughout planning and replaced by reconciliation.
# Retain those input records in the copy memo instead of cloning payloads
# that the detached finish pass immediately rebuilds. Geometry remains owned
# by the detached model and the callback keeps its original userdata.
function _reclassify_model_copy(model)
    memory=IdDict{Any,Any}()
    for record in values(model.discrete)
        memory[record]=record
    end
    for record in values(model.meshing.attached)
        memory[record]=record
    end
    return Model._model_mesh_transaction_copy(model,memory)
end

function _reclassify_row_maps(source,target,plan)
    node_positions=Dict(tag=>Int32(row) for (row,tag) in
        enumerate(target.entity_data.external_node_tags))
    target_cells=reduce(vcat,target.entity_data.external_element_tags;init=UInt64[])
    cell_positions=Dict(tag=>Int32(row) for (row,tag) in enumerate(target_cells))
    node_map=Int32[get(node_positions,get(plan.node_tag_map,tag,UInt64(0)),Int32(0))
        for tag in source.entity_data.external_node_tags]
    source_cells=reduce(vcat,source.entity_data.external_element_tags;init=UInt64[])
    cell_map=Int32[get(cell_positions,get(plan.element_tag_map,tag,UInt64(0)),Int32(0))
        for tag in source_cells]
    return node_map,cell_map
end

function _reclassify_remap_sections(source,node_map,cell_map,caller)
    sections=Elements.MshDataSection[]
    for section in source.data_sections
        nodes=Int32[];cells=Int32[];widths=Int32[];values=Float64[]
        node_offset=0;value_offset=0;components=Int(section.header[2])
        count=length(section.name=="NodeData" ? section.nodes : section.elements)
        for row in 1:count
            width=section.name=="ElementNodeData" ? Int(section.row_nodes[row]) : 1
            value_count=width*components
            selected=section.name=="NodeData" ? node_map[section.nodes[row]] :
                cell_map[section.elements[row]]
            if selected!=0
                if section.name=="NodeData"
                    push!(nodes,selected)
                else
                    push!(cells,selected)
                    if section.name=="ElementNodeData"
                        push!(widths,Int32(width))
                        if !section.implicit_nodes
                            for position in node_offset+1:node_offset+width
                                mapped=node_map[section.nodes[position]]
                                mapped!=0 || throw(ArgumentError(
                                    "$caller: retained ElementNodeData references a removed node"))
                                push!(nodes,mapped)
                            end
                        end
                    end
                end
                for position in value_offset+1:value_offset+value_count
                    push!(values,section.values[position])
                end
            end
            node_offset+=section.name=="ElementNodeData" ? width : 0
            value_offset+=value_count
        end
        header=copy(section.header)
        header[3]=length(section.name=="NodeData" ? nodes : cells)
        push!(sections,Elements.MshDataSection(section.name,copy(section.strings),
            copy(section.reals),header,cells,nodes,widths,section.implicit_nodes,
            values,section.anchor))
    end
    return sections
end

function _reclassify_data_sections(source,target,plan,caller)
    isempty(source.data_sections) && return Elements.MshDataSection[]
    node_map,cell_map=_reclassify_row_maps(source,target,plan)
    return _reclassify_remap_sections(source,node_map,cell_map,caller)
end

function _reclassify_partition_data(source,target,plan,caller)
    previous=source.partition_data
    previous===nothing && return nothing
    _,cell_map=_reclassify_row_maps(source,target,plan)
    ghosts=Elements.MixedGhostElement[]
    for ghost in previous.ghost_elements
        mapped=cell_map[ghost.element]
        mapped==0 && continue
        push!(ghosts,Elements.MixedGhostElement(mapped,ghost.partition,ghost.ghost_partitions))
    end
    return Elements.MixedPartitionData(previous.num_partitions;
        entities=previous.entities,ghost_entities=previous.ghost_entities,ghost_elements=ghosts)
end

# GModel::pruneMeshVertexAssociations associates every actual element vertex,
# including high-order supports, with the first entity at its lowest incident
# dimension. Cartesian positions do not participate in classification.
function _mesh_update_node_parameters(previous,supplied,caller)
    values=Float64[Float64(value) for value in supplied]
    all(isfinite,values) || throw(ArgumentError("$caller: parametric coordinates must be finite"))
    previous===nothing && return nothing
    parameters=copy(previous)
    for index in 1:min(length(parameters),length(values))
        parameters[index]=values[index]
    end
    return parameters
end

function _reclassify_element_nodes!(assembly;chains=Set{Tuple{Int,Int}}(),
                                    preserve_mesh_owners::Bool=false)
    owners=fill((4,Int32(0)),length(assembly.nodes))
    order=sortperm(assembly.cells;by=cell->(-cell.owner[1],cell.owner[2]),
        alg=Base.Sort.MergeSort)
    for position in order
        cell=assembly.cells[position]
        for node in cell.nodes
            owners[node][1]>cell.owner[1] && (owners[node]=cell.owner)
        end
    end
    # GModel first associates ordinary elements, then forces chain Regions,
    # Faces, Edges and Vertices in that order, with ascending entity tags.
    # The final chain owner can therefore override a lower ordinary owner.
    chain_order=sort!([position for position in eachindex(assembly.cells)
        if (assembly.cells[position].owner[1],Int(assembly.cells[position].owner[2])) in chains];
        by=position->(-assembly.cells[position].owner[1],assembly.cells[position].owner[2]),
        alg=Base.Sort.MergeSort)
    for position in chain_order
        cell=assembly.cells[position]
        for node in cell.nodes
            owners[node]=cell.owner
        end
    end
    orphan=falses(length(assembly.nodes))
    for position in eachindex(assembly.nodes)
        node=assembly.nodes[position]
        orphan[position]=owners[position][1]==4
        orphan[position] && continue
        preserve_mesh_owners ||
            (assembly.nodes[position]=_Dim01Node(node.tag,node.point,owners[position],node.parameter))
    end
    return orphan
end

function _reclassify_nodes_plan(model,cached,class,caller;node_max,element_max,
                                preserve_mesh_owners::Bool=false,
                                cell_edit=nothing)
    staged=_reclassify_model_copy(model)
    actual,actual_class=_dim01_actual_cache(cached,class;physical_names=model.physical_names)
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        node_max,element_max,typemax(Int32),typemax(Int32),caller)
    records=Model._discrete_mesh_records_model(staged)
    raw_nodes=sum(length(record.node_tags) for (_,record) in records;init=0)
    raw_cells=sum(length(record.element_tags) for (_,record) in records;init=0)
    node_capacity=raw_nodes+(actual===nothing ? 0 : nnodes(actual))
    cell_capacity=raw_cells+(actual===nothing ? 0 : _mesh_element_offsets(actual)[3])
    sizehint!(assembly.nodes,node_capacity)
    sizehint!(assembly.node_positions,node_capacity)
    sizehint!(assembly.cells,cell_capacity)
    # This operation changes associations, not interpolation order. Retain
    # every admitted source MSH slot, including higher-order supports.
    source_rows=_dim01_ingest!(assembly,staged,actual,actual_class;
                              preserve_source_order=true)
    # Surgical plugins (Crack) rewrite connectivity after ingest: `cell_edit`
    # receives the assembly and may retarget `cell.nodes` positions, e.g. to
    # duplicated nodes staged on a fresh record.
    cell_edit===nothing || cell_edit(assembly)
    orphan=_reclassify_element_nodes!(assembly;chains=staged.meshing.chain_entities,
                                      preserve_mesh_owners)
    sort!(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
    plan=_dim01_finish!(assembly,staged,model_entities(staged),
        Set{Tuple{Int,Int32}}(),orphan,false,false;
        source_mesh=actual,source_node_positions=source_rows,node_tag_order=true)
    mesh=plan.mesh
    if actual isa MixedMesh && (!isempty(actual.ancillary_sections) ||
            !isempty(actual.data_sections) || actual.partition_data!==nothing ||
            actual.physical_names!=mesh.physical_names)
        # Preserve imported names while including new physical groups from
        # the detached model (e.g. freshly stored homology chains).
        physical_names=merge(actual.physical_names,mesh.physical_names)
        preserved=MixedMesh(mesh.coords,mesh.blocks;physical_names,
            entity_data=mesh.entity_data,elementary_entities=mesh.elementary_entities,
            periodic_links=mesh.periodic_links,ancillary_sections=actual.ancillary_sections,
            data_sections=_reclassify_data_sections(actual,mesh,plan,caller),
            partition_data=_reclassify_partition_data(actual,mesh,plan,caller))
        plan=merge(plan,(mesh=preserved,class=_mixed_rebind_class(plan.class,preserved),))
    end
    previous=staged.meshing.mesh_authority
    previous===nothing || _reclassify_rebind_carrier!(plan,previous,caller)
    return plan
end

function _reclassify_nodes()
    caller="API.mesh.reclassify_nodes"
    return lock(STATE_LOCK) do
        model=_model_locked()
        cached=LAST_MESH[]
        class=cached===nothing ? nothing : _cached_classification_locked(cached)
        plan=_reclassify_nodes_plan(model,cached,class,caller;
            node_max=NODE_TAG_MAX[],element_max=ELEMENT_TAG_MAX[])
        # _dim01_finish! has already reconciled these exact public labels
        # and original vertex parameters on the detached model.
        _commit_tagged_plan_locked!(plan;preserve_visibility=true,reconcile_records=false)
        return nothing
    end
end

function _reclassify_rebind_carrier!(plan,previous,caller)
    carrier=previous.carrier
    (carrier===nothing || carrier===previous.mesh) &&
        return _higher_rebind_plan_authority!(plan,previous)
    positions=Dict(tag=>Int32(row) for (row,tag) in
        enumerate(plan.mesh.entity_data.external_node_tags))
    rows=Int32[get(positions,get(plan.node_tag_map,
        previous.mesh.entity_data.external_node_tags[old],UInt64(0)),Int32(0))
        for old in previous.carrier_rows]
    kept=findall(!iszero,rows)
    mapping=zeros(Int32,length(rows));mapping[kept]=Int32.(1:length(kept))
    if carrier isa Mesh
        all(node->mapping[node]!=0,carrier.segs) &&
            all(node->mapping[node]!=0,carrier.tris) &&
            all(node->mapping[node]!=0,carrier.tets) || throw(ArgumentError(
                "$caller: native carrier connectivity has no retained actual identity"))
        retained=Mesh(carrier.coords[:,kept];segs=mapping[carrier.segs],
            tris=mapping[carrier.tris],tets=mapping[carrier.tets],
            seg_tag=carrier.seg_tag,tri_tag=carrier.tri_tag,tet_tag=carrier.tet_tag)
    else
        blocks=ElementBlock[]
        for block in carrier.blocks
            all(node->mapping[node]!=0,block.nodes) || throw(ArgumentError(
                "$caller: native carrier connectivity has no retained actual identity"))
            push!(blocks,ElementBlock(block.msh,mapping[block.nodes],block.tags))
        end
        previous_data=carrier.entity_data
        coordinates=copy(plan.mesh.coords[:,rows[kept]])
        parameters=deepcopy(plan.mesh.entity_data.node_parametric[rows[kept]])
        for (slot,old_slot) in enumerate(kept)
            old_row=previous.carrier_rows[old_slot];new_row=rows[old_slot]
            unchanged=previous.mesh.entity_data.node_entities[old_row]==plan.mesh.entity_data.node_entities[new_row] &&
                isequal(previous.mesh.coords[:,old_row],plan.mesh.coords[:,new_row]) &&
                isequal(previous.mesh.entity_data.node_parametric[old_row],
                        plan.mesh.entity_data.node_parametric[new_row])
            if unchanged
                coordinates[:,slot]=carrier.coords[:,old_slot]
                previous_data===nothing || (parameters[slot]=previous_data.node_parametric[old_slot])
            end
        end
        data=previous_data===nothing ? nothing : Elements._native_mixed_entity_data(
            plan.mesh.entity_data.entities;node_entities=plan.mesh.entity_data.node_entities[rows[kept]],
            node_parametric=parameters,
            external_node_tags=plan.mesh.entity_data.external_node_tags[rows[kept]],
            block_entities=previous_data.block_entities,
            external_element_tags=previous_data.external_element_tags)
        cells=sum(size(block.nodes,2) for block in carrier.blocks;init=0)
        retained=MixedMesh(coordinates,blocks;
            physical_names=carrier.physical_names,entity_data=data,
            elementary_entities=carrier.elementary_entities,
            periodic_links=_mixed_remap_periodic_nodes(carrier,mapping,caller;drop_removed=true),
            ancillary_sections=carrier.ancillary_sections,
            data_sections=_reclassify_remap_sections(carrier,mapping,collect(Int32,1:cells),caller),
            partition_data=carrier.partition_data)
    end
    cells=Model._model_mesh_rebound_cell_tags(previous,plan.element_tag_map)
    return _higher_rebind_plan_authority!(plan,previous;
        carrier=retained,carrier_rows=rows[kept],carrier_cells=cells)
end
