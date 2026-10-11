# A homology heuristic uses the actual MElement nodal map. Resolve cache,
# native overlay and record vertices once in the same public label space as
# its cell connectivity, including records with cache-only references.
function _homology_public_coordinates(model,cached,class)
    coordinates=Dict{UInt64,NTuple{3,Float64}}()
    for (_,record) in Model._discrete_mesh_records_model(model)
        for (row,tag) in enumerate(record.node_tags)
            coordinates[tag]=ntuple(axis->record.node_coords[axis,row],3)
        end
    end
    if cached!==nothing
        overlay=_high_order_overlay(cached)
        count=overlay===nothing ? nnodes(cached) : nnodes(overlay)
        public=class===nothing ? nothing : class.public_tags
        tags=_cache_node_tags(public,1:count)
        for (row,tag) in enumerate(tags)
            geometry=row<=nnodes(cached) ? cached : overlay
            coordinates[UInt64(tag)]=ntuple(axis->geometry.coords[axis,row],3)
        end
    end
    return coordinates
end
