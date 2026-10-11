function _geo_restore_orientation_constraints!(m)
    for (tag,flag) in m.meshing.orientation_geo_reverse
        haskey(m.surfaces,tag) || continue
        m.meshing.reverse[(2,tag)]=flag
    end
    return nothing
end

function _geo_orientation_source_flags!(m,tags)
    for region in tags,(dimension,tag) in Model.model_boundary(m,[(3,region)],false,false,false)
        dimension==2 || continue
        haskey(m.discrete,(2,tag)) && continue
        Model._model_orientation_occ_plane(m,tag) && continue
        geometry=get(m.surface_geometry,tag,nothing)
        geometry!==nothing && hasproperty(geometry,:occ) && continue
        get!(m.meshing.orientation_geo_reverse,tag,
             get(m.meshing.reverse,(2,tag),false))
    end
    return nothing
end

function _geo_orientation_blocks(mesh::Mesh)
    return ((1,mesh.segs,mesh.seg_tag),(2,mesh.tris,mesh.tri_tag),
            (4,mesh.tets,mesh.tet_tag))
end
function _geo_orientation_blocks(mesh::MixedMesh)
    return [(Int(block.msh),block.nodes,block.tags) for block in mesh.blocks]
end

# Legacy merge columns follow explicit part/family order or mixed bucket order.
# Stored node remaps identify geometry rows even after coordinate edits.
function _geo_orientation_legacy_input(m,context,regions)
    surfaces,volumes=Model._model_orientation_entities(m,regions)
    mesh=context.mesh
    mesh===nothing && return Model._model_orientation_record_input(m;
        surface_tags=surfaces,volume_tags=volumes),
        Dict{Int,Tuple{Int,Int,Int,Int}}()
    _geo_check_part_rows(context)
    length(context.mesh_node_parametric)==nnodes(mesh) || throw(ArgumentError(
        "execute_geo: ReorientMesh stored parameters do not align with mesh rows"))
    owners=[get(context.mesh_node_owner,row,(-1,0)) for row in 1:nnodes(mesh)]
    count=sum(size(cells,2) for (_,cells,_) in _geo_orientation_blocks(mesh);init=0)
    input=Model._model_orientation_record_input(m;coords=mesh.coords,
        node_tags=UInt64.(1:nnodes(mesh)),owners,parameters=context.mesh_node_parametric,
        key_offset=count,surface_tags=surfaces,volume_tags=volumes)
    bindings=Dict{Int,Tuple{Int,Int,Int,Int}}()
    bucket_indices=Dict{Tuple{Int,Int,Int},Int}()
    bucket_columns=Dict{Tuple{Int,Int,Int},Int}()
    family_columns=Dict(1=>0,2=>0,4=>0)
    block_offsets=Int[];offset=0
    for (_,cells,_) in _geo_orientation_blocks(mesh)
        push!(block_offsets,offset);offset+=size(cells,2)
    end
    for (part_index,((dimension,tag,part),rows)) in
            enumerate(zip(context.mesh_parts,context.mesh_part_rows))
        for (local_block,(msh,cells,_)) in enumerate(_geo_orientation_blocks(part))
            isempty(cells) && continue
            bucket=(msh,dimension,tag)
            block=mesh isa Mesh ? msh==1 ? 1 : msh==2 ? 2 : 3 :
                get!(bucket_indices,bucket,length(bucket_indices)+1)
            base=mesh isa Mesh ? family_columns[msh] : get(bucket_columns,bucket,0)
            for column in axes(cells,2)
                key=block_offsets[block]+base+column
                if dimension==2 && msh_dimension(msh)==2 && tag in surfaces
                    push!(get!(()->Model._ModelMeshOrientationCell[],input.surface_cells,tag),
                        Model._ModelMeshOrientationCell(key,msh,rows[cells[:,column]]))
                    bindings[key]=(part_index,local_block,column,base+column)
                elseif dimension==3 && msh_dimension(msh)==3 && tag in volumes
                    append!(get!(()->Int32[],input.volume_nodes,tag),rows[cells[:,column]])
                end
            end
            if mesh isa Mesh
                family_columns[msh]+=size(cells,2)
            else
                bucket_columns[bucket]=base+size(cells,2)
            end
        end
    end
    return input,bindings
end

function _geo_orientation_replace(mesh,changes)
    mesh===nothing && return nothing
    mesh isa MixedMesh && return Model._model_orientation_mixed_replacement(mesh,changes)
    triangles=copy(mesh.tris);changed=false
    for column in axes(triangles,2)
        nodes=get(changes,nsegs(mesh)+column,nothing)
        nodes===nothing && continue
        nodes==@view(triangles[:,column]) && continue
        triangles[:,column]=nodes;changed=true
    end
    changed || return mesh
    return Mesh(mesh.coords;segs=mesh.segs,tris=triangles,tets=mesh.tets,
        seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,tet_tag=mesh.tet_tag)
end

function _geo_orientation_record_changes!(input,plan)
    inverse=Dict(row=>Int32(node) for (node,row) in input.rows)
    for (key,(record,index)) in input.bindings
        nodes=get(plan.cells,key,nothing)
        nodes===nothing || (record.element_nodes[index]=Int32[inverse[row] for row in nodes])
    end
    return nothing
end

function _geo_orientation_actual!(m,context,regions,caller)
    working=_geo_mesh_working_model(m,context)
    _geo_orientation_source_flags!(working,regions)
    previous=working.meshing.mesh_authority
    actual=previous.mesh
    covered=Set((dimension,Int(tag)) for (dimension,tag) in keys(previous.products))
    input=Model._model_orientation_mixed_input(working,actual;covered,regions)
    plan=Model._model_prepare_mesh_orientation(working,regions,input.coords,
        input.surface_cells,input.owners,input.parameters;volume_nodes=input.volume_nodes,
        node_tags=input.node_tags,
        old_ruled_surface=_geo_exec_old_ruled(context))
    _geo_orientation_record_changes!(input,plan)
    replacement=Model._model_orientation_mixed_replacement(actual,plan.cells)
    carrier=Model._model_orientation_carrier(actual,replacement,previous.carrier,
                                            previous.carrier_rows,plan.cells;cell_tags=previous.carrier_cells)
    working.meshing.mesh_authority=Model._model_orientation_authority(previous,replacement,carrier)
    Model._model_publish_mesh_orientation!(working,plan)
    completed=(;model=working,mesh=replacement,max_node_tag=context.mesh_node_max,
               max_element_tag=context.mesh_element_max)
    _geo_publish_actual_plan!(m,context,completed,caller;carrier)
    for message in plan.diagnostics;_geo_yywarn!(context,message);end
    return nothing
end

function _geo_exec_reorient_mesh!(m,list_source,context)
    caller="execute_geo: ReorientMesh Volume"
    list=_geo_constraint_list(list_source,context,caller;allow_all=false)
    _geo_sync_physical_view_if_changed!(m,context)
    regions=Int[]
    for value in list
        tag=_geo_constraint_int(value,caller,"Volume tag")
        _has_entity(m,3,tag) && push!(regions,tag)
    end
    isempty(regions) && return nothing
    if hasproperty(m.meshing,:mesh_authority)
        authority=getproperty(m.meshing,:mesh_authority)
        if authority!==nothing && authority.mesh!==nothing
            return _geo_orientation_actual!(m,context,regions,caller)
        end
    end
    working=Model._model_mesh_transaction_copy(m)
    _geo_orientation_source_flags!(working,regions)
    input,bindings=_geo_orientation_legacy_input(working,context,regions)
    plan=Model._model_prepare_mesh_orientation(working,regions,input.coords,
        input.surface_cells,input.owners,input.parameters;volume_nodes=input.volume_nodes,
        node_tags=input.node_tags,
        old_ruled_surface=_geo_exec_old_ruled(context))
    _geo_orientation_record_changes!(input,plan)
    replacement=_geo_orientation_replace(context.mesh,plan.cells)
    parts=copy(context.mesh_parts)
    prepared=Dict{Int,Union{Mesh,MixedMesh}}()
    for (key,(part_index,local_block,column,_)) in bindings
        nodes=get(plan.cells,key,nothing)
        nodes===nothing && continue
        part=get!(prepared,part_index) do
            deepcopy(parts[part_index][3])
        end
        rows=context.mesh_part_rows[part_index]
        inverse=Dict(row=>Int32(index) for (index,row) in enumerate(rows))
        _,cells,_=_geo_orientation_blocks(part)[local_block]
        cells[:,column]=Int32[inverse[row] for row in nodes]
    end
    for (index,part) in prepared
        dimension,tag,_=parts[index];parts[index]=(dimension,tag,part)
    end
    Model._model_publish_mesh_orientation!(working,plan)
    for name in fieldnames(GeoModel);setfield!(m,name,getfield(working,name));end
    context.mesh=replacement;context.mesh_parts=parts
    for message in plan.diagnostics;_geo_yywarn!(context,message);end
    return nothing
end
