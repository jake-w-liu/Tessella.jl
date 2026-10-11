# Imported Ghost collections contain the elements assigned by GhostElements;
# the full imported element objects remain separately available for UID lookup.
function _public_bulk_ghost_positions(model)
    authority=model.meshing.mesh_authority
    mesh=authority===nothing ? nothing : authority.mesh
    positions=Dict{UInt64,Int}()
    if mesh!==nothing && mesh.entity_data!==nothing
        offset=0
        for tags in mesh.entity_data.external_element_tags
            for uid in tags
                offset+=1;positions[uid]=offset
            end
        end
    end
    result=Dict{Tuple{Int,Int},Vector{Int}}()
    for (key,receipt) in model.meshing.imported_entities
        startswith(receipt.kind,"Ghost ") || continue
        members=Int[]
        for (uid,_) in get(model.meshing.imported_ghost_elements,key,Tuple{UInt64,Int32}[])
            position=get(positions,uid,0)
            iszero(position) || push!(members,position)
        end
        result[key]=members
    end
    return result
end

@inline function _public_bulk_entity_visible(model,dimension,tag;ghost_positions=nothing)
    receipt=get(model.meshing.imported_entities,(Int(dimension),Int(tag)),nothing)
    (receipt===nothing || !startswith(receipt.kind,"Ghost ")) && return true
    ghosts=ghost_positions===nothing ? _public_bulk_ghost_positions(model) : ghost_positions
    return !isempty(get(ghosts,(Int(dimension),Int(tag)),Int[]))
end

function _native_has_ghost_entities(model)
    return any(receipt->startswith(receipt.kind,"Ghost "),
        values(model.meshing.imported_entities))
end

function _imported_ghost_uid_owners(model)
    owners=Dict{UInt64,Vector{Tuple{Int,Int32}}}()
    for key in sort!(collect(keys(model.meshing.imported_ghost_elements)))
        receipt=get(model.meshing.imported_entities,key,nothing)
        receipt!==nothing && startswith(receipt.kind,"Ghost ") || continue
        for (uid,_) in model.meshing.imported_ghost_elements[key]
            push!(get!(()->Tuple{Int,Int32}[],owners,uid),(key[1],Int32(key[2])))
        end
    end
    return owners
end

function _imported_ghost_uid_ranks(model)
    return Dict(key=>Dict(uid=>position for (position,(uid,_)) in enumerate(members))
        for (key,members) in model.meshing.imported_ghost_elements)
end

function _public_bulk_cell_positions(model,dimension,owners,entity;offset=0,ghost_positions=nothing)
    ghost_positions===nothing && !_native_has_ghost_entities(model) &&
        return findall(owner->entity<0 || owner==entity,owners)
    ghosts=ghost_positions===nothing ? _public_bulk_ghost_positions(model) : ghost_positions
    selected=entity<0 ? sort!(unique!(vcat(Int.(owners),
        [tag for ((dim,tag),_) in ghosts if dim==dimension]))) : Int[entity]
    result=Int[]
    for tag in selected
        key=(Int(dimension),tag)
        if haskey(ghosts,key)
            for position in ghosts[key]
                offset<position<=offset+length(owners) && push!(result,position-offset)
            end
        else
            append!(result,findall(==(tag),owners))
        end
    end
    return result
end

function _public_bulk_family_positions(model,class,msh,entity)
    result=Int[];family_offset=0
    ghost_positions=_native_has_ghost_entities(model) ? _public_bulk_ghost_positions(model) : nothing
    for (kind,dimension,offset,cells,owners) in _cache_catalog(class.mesh,class)
        kind==msh || continue
        positions=_public_bulk_cell_positions(model,dimension,owners,entity;offset,ghost_positions)
        append!(result,family_offset .+ positions)
        family_offset+=size(cells,2)
    end
    return result
end

function _imported_unpartition_plan(model,cached,class,caller)
    receipts=model.meshing.imported_entities
    removed=Set{Tuple{Int,Int}}(key for (key,value) in receipts
        if startswith(value.kind,"Partition ") || startswith(value.kind,"Ghost "))
    isempty(removed) && return nothing
    cached isa MixedMesh && cached.entity_data!==nothing && class!==nothing ||
        throw(ArgumentError("$caller: imported partition entities require their actual mesh authority"))
    assembly=_Dim01Assembly(_Dim01Node[],_Dim01Cell[],Dict{UInt64,Int32}(),
        NODE_TAG_MAX[],ELEMENT_TAG_MAX[],typemax(Int32),typemax(Int32),caller)
    _dim01_ingest!(assembly,model,cached,class;preserve_source_order=true)
    visited=falses(length(assembly.nodes))
    transferred=Int[]
    retained=_Dim01Cell[]
    for position in sortperm(assembly.cells;by=_dim01_cell_order,alg=Base.Sort.MergeSort)
        cell=assembly.cells[position]
        key=(cell.owner[1],Int(cell.owner[2]))
        if key in removed
            receipt=receipts[key]
            startswith(receipt.kind,"Ghost ") && continue
            parent=receipt.parent
            (parent[1]==key[1] && Model._model_entity_known(model,parent...)) || continue
            parent in removed && throw(ArgumentError(
                "$caller: a partition parent is itself removed by unpartitioning"))
            cell.owner=(parent[1],Int32(parent[2]))
            for row in cell.nodes
                visited[row] && continue
                visited[row]=true
                node=assembly.nodes[row]
                assembly.nodes[row]=_Dim01Node(node.tag,node.point,cell.owner,node.parameter)
                push!(transferred,row)
            end
        end
        push!(retained,cell)
    end
    assembly.cells=retained
    # Parent vertices are appended in first element-incidence order, after
    # the vertices that the parent already owned. No Cartesian search occurs.
    order=vcat(findall(!,visited),transferred)
    remap=zeros(Int32,length(order));remap[order]=Int32.(1:length(order))
    assembly.nodes=assembly.nodes[order]
    empty!(assembly.node_positions)
    for (row,node) in enumerate(assembly.nodes)
        assembly.node_positions[node.tag]=Int32(row)
    end
    for cell in assembly.cells
        cell.nodes=remap[cell.nodes]
        any(row->(assembly.nodes[row].owner[1],Int(assembly.nodes[row].owner[2])) in removed,
            cell.nodes) && throw(ArgumentError(
                "$caller: a surviving element references an untransferred partition vertex"))
    end
    staged=Model._model_mesh_transaction_copy(model)
    state=Model._model_removal_state(staged,removed)
    for (name,value) in pairs(state)
        setfield!(staged,name,value)
    end
    staged.meshing.imported_num_partitions=0
    # Native partition-group names contain underscores. UnpartitionMesh
    # removes these groups while retaining every other surviving membership.
    for (key,name) in collect(staged.physical_names)
        occursin('_',name) || continue
        delete!(staged.physical,key);delete!(staged.physical_names,key)
        for (entity,groups) in staged.entity_physicals
            entity[1]==key[1] && filter!(!=(key[2]),groups)
        end
    end
    for cell in assembly.cells
        cell.physical=_dim01_physical(staged,(cell.owner[1],Int(cell.owner[2])))
    end
    discarded=Set{Tuple{Int,Int32}}((d,Int32(t)) for (d,t) in removed)
    plan=_dim01_finish!(assembly,staged,model_entities(staged),discarded,
        falses(length(assembly.nodes)),false,false;record_view=true)
    target=plan.mesh
    node_map,cell_map=_reclassify_row_maps(cached,target,plan)
    links=Elements.MixedPeriodicLink[]
    for link in cached.periodic_links
        ((link.dim,Int(link.slave_entity)) in removed ||
         (link.dim,Int(link.master_entity)) in removed) && continue
        slaves=node_map[link.slave_nodes];masters=node_map[link.master_nodes]
        keep=findall(i->slaves[i]!=0 && masters[i]!=0,eachindex(slaves))
        push!(links,Elements.MixedPeriodicLink(link.dim,link.slave_entity,link.master_entity,
            slaves[keep],masters[keep];affine=link.affine))
    end
    sections=_reclassify_remap_sections(cached,node_map,cell_map,caller)
    target=MixedMesh(target.coords,target.blocks;physical_names=target.physical_names,
        entity_data=target.entity_data,elementary_entities=target.elementary_entities,
        periodic_links=links,ancillary_sections=cached.ancillary_sections,data_sections=sections)
    previous=model.meshing.mesh_authority
    ranks=previous===nothing || previous.mesh===nothing ? Dict{UInt64,UInt64}() :
        Dict(tag=>rank for (tag,rank) in zip(previous.mesh.entity_data.external_node_tags,
                                           previous.private_rank))
    next=maximum(values(ranks);init=UInt64(0))
    private_rank=UInt64[]
    for tag in target.entity_data.external_node_tags
        rank=get(ranks,tag,UInt64(0))
        if iszero(rank)
            next<typemax(UInt64) || throw(ArgumentError("$caller: private vertex ranks overflow UInt64"))
            next+=UInt64(1);rank=next
        end
        push!(private_rank,rank)
    end
    staged.meshing.mesh_authority=Model._model_mesh_authority(target,plan.authority;
        caller,private_rank)
    rebound=_mixed_rebind_class(plan.class,target)
    return merge(plan,(mesh=target,class=rebound,preserved_overlay=nothing,
        preserved_mid_owners=Tuple{Int,Int32}[]))
end
