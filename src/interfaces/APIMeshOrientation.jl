# Outward orientation prepares metadata and connectivity before publication.
# Actual P2 support coordinates and public row identities are retained.
function _orientation_cache_input(model,cached,class,overlay,regions,caller)
    surfaces,volumes=Model._model_orientation_entities(model,regions)
    if cached===nothing
        return Model._model_orientation_record_input(model;
            surface_tags=surfaces,volume_tags=volumes)
    end
    class===nothing && throw(ArgumentError(
        "$caller: existing cache requires explicit entity classification"))
    coordinates=overlay===nothing ? cached.coords : overlay.coords
    public=class.public_tags
    tags=public===nothing ? UInt64.(1:size(coordinates,2)) : public.node_tags
    owners=Tuple{Int,Int}[(d,Int(t)) for (d,t) in class.node_entities]
    parameters=copy(class.node_parametric)
    if overlay!==nothing
        append!(owners,Tuple{Int,Int}[(d,Int(t)) for (d,t) in LAST_MESH_HIGH_ORDER_MIDS[]])
        for node in nnodes(cached)+1:nnodes(overlay)
            push!(parameters,get(class.mid_parametric,_relocation_mid_key(class,node),nothing))
        end
    end
    count=last(_mesh_element_offsets(cached))
    records=[((d,t),record) for (d,t,record) in _discrete_mesh_records(model)]
    input=Model._model_orientation_record_input(model;coords=coordinates,
        node_tags=tags,owners,parameters,key_offset=count,records,
        surface_tags=surfaces,volume_tags=volumes)
    record_owners=IdDict(record=>owner for (owner,record) in records)
    for (key,(record,index)) in input.bindings
        owner=record_owners[record]
        target=_raw_model_mesh_record(model,owner...)
        position=findfirst(==(record.element_tags[index]),target.element_tags)
        position===nothing && throw(ArgumentError(
            "$caller: retained element has no stored record identity"))
        input.bindings[key]=(target,position)
    end
    for (msh,dimension,offset,stored,cell_owners) in _cache_catalog(cached,class)
        dimension in (2,3) || continue
        elevated=overlay!==nothing && msh==_p2_skeleton_etype(overlay)
        cells=elevated ? _p2_cells(overlay) : stored
        kind=elevated ? Int(_p2_etype(overlay)) : Int(msh)
        for column in axes(cells,2)
            tag=Int(cell_owners[column])
            tag in (dimension==2 ? surfaces : volumes) || continue
            if dimension==2
                push!(get!(()->Model._ModelMeshOrientationCell[],input.surface_cells,tag),
                    Model._ModelMeshOrientationCell(offset+column,kind,
                        @view(cells[:,column])))
            else
                append!(get!(()->Int32[],input.volume_nodes,tag),@view(cells[:,column]))
            end
        end
    end
    return input
end

function _orientation_replacement(cached,class,overlay,changes,caller)
    cached===nothing && return nothing,nothing,nothing
    replacement=cached
    changed=false
    if cached isa MixedMesh
        blocks=ElementBlock[]
        offset=0
        for block in cached.blocks
            cells=copy(block.nodes)
            for column in axes(cells,2)
                nodes=get(changes,offset+column,nothing)
                nodes===nothing && continue
                nodes==@view(cells[:,column]) && continue
                cells[:,column]=nodes;changed=true
            end
            push!(blocks,ElementBlock(block.msh,cells,block.tags))
            offset+=size(cells,2)
        end
        changed && (replacement=_mixed_with_coordinates(cached,cached.coords;blocks))
        return replacement,changed ? _mixed_rebind_class(class,replacement) : class,nothing
    end
    triangles=copy(cached.tris)
    offset=nsegs(cached)
    for column in axes(triangles,2)
        nodes=get(changes,offset+column,nothing)
        nodes===nothing && continue
        @view(nodes[1:3])==@view(triangles[:,column]) && continue
        triangles[:,column]=nodes[1:3];changed=true
    end
    changed && (replacement=Mesh(cached.coords;segs=cached.segs,tris=triangles,
        tets=cached.tets,seg_tag=cached.seg_tag,tri_tag=cached.tri_tag,tet_tag=cached.tet_tag))
    public=class.public_tags
    table=!changed || public===nothing ? public : _cache_public_tags(replacement;
        authority=public.authority,node_tags=public.node_tags,element_tags=public.element_tags,
        node_count=length(public.node_tags),caller)
    rebound=changed ? _MeshClassification(replacement,class.entity,class.entities,
        class.node_entities,class.boundaries,class.seg_entities,class.tri_entities,
        class.tet_entities,class.cell_entities,table,class.edge_entities,class.face_entities,
        class.quad_entities,class.node_parametric,copy(class.mid_parametric),copy(class.mid_keys)) : class
    prepared=nothing
    if overlay!==nothing
        actual=overlay
        if changed && overlay isa P2TriMesh
            cells=copy(_p2_cells(overlay))
            for column in axes(cells,2)
                nodes=get(changes,offset+column,nothing)
                nodes===nothing || (cells[:,column]=nodes)
            end
            actual=P2TriMesh(overlay.coords,cells;tri_tag=overlay.tri_tag,_check_geometry=false)
        end
        prepared=(actual,copy(LAST_MESH_HIGH_ORDER_MIDS[]))
    end
    return replacement,rebound,prepared
end

function _orient_mesh_outward(tag)
    caller="API.mesh.set_outward_orientation"
    return lock(STATE_LOCK) do
        original=_model_locked()
        region=Model._tag(tag,caller,3)
        Model._mesh_attr_entity!(original,3,region,caller;allow_discrete=true)
        staged=Model._model_mesh_transaction_copy(original)
        cached=LAST_MESH[]
        class=cached===nothing ? nothing : _cached_classification_locked(cached)
        overlay=_high_order_overlay(cached)
        input=_orientation_cache_input(staged,cached,class,overlay,[region],caller)
        plan=Model._model_prepare_mesh_orientation(staged,[region],input.coords,
            input.surface_cells,input.owners,input.parameters;volume_nodes=input.volume_nodes,
            node_tags=input.node_tags,
            old_ruled_surface=!iszero(OPTIONS["Geometry.OldRuledSurface"]))
        inverse=Dict(row=>UInt64(node) for (node,row) in input.rows)
        for (key,(record,index)) in input.bindings
            nodes=get(plan.cells,key,nothing)
            nodes===nothing || (record.element_nodes[index]=UInt64[inverse[row] for row in nodes])
        end
        replacement,rebound,prepared=_orientation_replacement(cached,class,overlay,
                                                               plan.cells,caller)
        mid_parameters=prepared===nothing ? nothing :
            _relocation_prepare_mid_parameters(staged,replacement,rebound,prepared[1],prepared[2])
        partition=LAST_MESH_PARTITION[]
        partition=partition===nothing || partition.mesh!==cached ? nothing :
            _MeshPartition(replacement,partition.num_partitions,partition.element_partitions)
        if hasproperty(original.meshing,:mesh_authority)
            previous=getproperty(original.meshing,:mesh_authority)
            if previous!==nothing && previous.mesh===cached
                carrier=Model._model_orientation_carrier(cached,replacement,previous.carrier,
                    previous.carrier_rows,plan.cells;cell_tags=previous.carrier_cells)
                staged.meshing.mesh_authority=Model._model_orientation_authority(
                    previous,replacement,carrier)
            elseif previous!==nothing && previous.carrier===cached && previous.mesh!==nothing
                throw(ArgumentError(
                    "$caller: a legacy geometry carrier requires the full actual-cache orientation adapter"))
            end
        end
        Model._model_publish_mesh_orientation!(staged,plan)
        # Geometry evaluation, record resolution, cache construction and stored
        # parameter preparation all finish before the following assignments.
        index=_current_slot_index_locked()
        index==0 || (MODEL_SLOTS[index].model=staged)
        CURRENT[]=staged
        LAST_MESH[]=replacement
        LAST_MESH_CLASS[]=rebound
        LAST_MESH_PARTITION[]=partition
        if replacement!==cached
            LAST_MESH_LOCATOR[]=nothing
            LAST_MIXED_MESH_LOCATOR[]=nothing
            LAST_MESH_EDGES[]=nothing
            LAST_MESH_FACES[]=nothing
        end
        prepared===nothing || _set_high_order_overlay(prepared[1],replacement,
            prepared[2];mid_parameters)
        for message in plan.diagnostics;@warn message;end
        return nothing
    end
end
