# Public labels admit a current owned product; coordinates only verify that
# explicit binding. Naked meshes retain the geometric projection contract.
function _model_mesh_projection_input(m::GeoModel,mesh,dimension::Int,tag::Int,caller)
    mesh isa MixedMesh && mesh.entity_data!==nothing || return nothing
    authority=m.meshing.mesh_authority
    authority!==nothing && authority.mesh!==nothing || return nothing
    owner=(dimension,Int32(tag))
    product=get(authority.products,owner,nothing)
    product!==nothing || throw(ArgumentError(
        "$caller: input does not identify a current actual entity product"))
    actual=authority.mesh;stored=actual.entity_data;data=mesh.entity_data
    positions=Dict(stored.external_node_tags[row]=>row for row in product.rows)
    length(positions)==nnodes(mesh) || throw(ArgumentError(
        "$caller: input node set does not match the current actual entity product"))
    rows=Vector{Int32}(undef,nnodes(mesh))
    for row in eachindex(rows)
        bound=get(positions,data.external_node_tags[row],Int32(0))
        bound>0 || throw(ArgumentError(
            "$caller: input node label is absent from the current actual entity product"))
        mesh.coords[:,row]==actual.coords[:,bound] || throw(ArgumentError(
            "$caller: input coordinates do not match their current actual node labels"))
        data.node_entities[row]==stored.node_entities[bound] &&
            data.node_parametric[row]==stored.node_parametric[bound] || throw(ArgumentError(
            "$caller: input ownership or parameters do not match their current actual node labels"))
        rows[row]=bound
        node_owner=data.node_entities[row]
        supplied=get(data.entities,node_owner,nothing)
        current=get(stored.entities,node_owner,nothing)
        supplied!==nothing && current!==nothing || throw(ArgumentError(
            "$caller: input node owner lacks its current actual entity metadata"))
        for field in fieldnames(MixedEntity)
            getfield(supplied,field)==getfield(current,field) || throw(ArgumentError(
                "$caller: input entity metadata does not match its current actual owner"))
        end
    end
    supplied=get(data.entities,owner,nothing);current=get(stored.entities,owner,nothing)
    supplied!==nothing && current!==nothing || throw(ArgumentError(
        "$caller: input cell owner lacks its current actual entity metadata"))
    for field in fieldnames(MixedEntity)
        getfield(supplied,field)==getfield(current,field) || throw(ArgumentError(
            "$caller: input entity metadata does not match its current actual owner"))
    end
    slots=Dict{UInt64,Tuple{Int,Int}}()
    for (block_index,tags) in enumerate(product.element_tags),column in eachindex(tags)
        slots[tags[column]]=(block_index,column)
    end
    seen=Set{UInt64}()
    for (block_index,block) in enumerate(mesh.blocks)
        block isa ElementBlock && msh_dimension(block.msh)==dimension || throw(ArgumentError(
            "$caller: actual input must contain only the selected entity's cell dimension"))
        for column in axes(block.nodes,2)
            data.block_entities[block_index][column]==tag || throw(ArgumentError(
                "$caller: input cell owner does not match the selected actual entity"))
            label=data.external_element_tags[block_index][column]
            slot=get(slots,label,nothing)
            slot!==nothing || throw(ArgumentError(
                "$caller: input cell label is absent from the current actual entity product"))
            reference=product.blocks[slot[1]];reference_column=slot[2]
            block.msh==reference.msh && block.tags[column]==reference.tags[reference_column] &&
                size(block.nodes,1)==size(reference.nodes,1) || throw(ArgumentError(
                "$caller: input cell family or physical tag does not match its current actual label"))
            for corner in axes(block.nodes,1)
                rows[block.nodes[corner,column]]==
                    product.rows[reference.nodes[corner,reference_column]] || throw(ArgumentError(
                    "$caller: input cell topology or winding does not match its current actual label"))
            end
            push!(seen,label)
        end
    end
    length(seen)==length(slots)>0 || throw(ArgumentError(
        "$caller: input cells do not match the complete current actual entity product"))
    return (authority=authority,rows=rows)
end

function _model_mesh_projection_closure(m::GeoModel,dimension::Int,tag::Int,caller)
    closure=Set{Tuple{Int,Int32}}(((dimension,Int32(tag)),))
    if dimension==3
        points,curves,surfaces,_,_,_=_model_volume_embedding_inventory(m,tag,caller)
        union!(closure,((2,Int32(surface)) for surface in surfaces))
        union!(closure,((1,Int32(curve)) for curve in curves))
        union!(closure,((0,Int32(point)) for point in points))
    else
        points,embedded=_model_surface_embedding_tags(m,tag,caller)
        curves=vcat(_model_projection_surface_curves(m,tag),embedded)
        for curve in curves
            push!(closure,(1,Int32(curve)))
            for point in m.curves[curve]
                push!(closure,(0,Int32(point)))
            end
        end
        union!(closure,((0,Int32(point)) for point in points))
    end
    return closure
end

function _model_mesh_actual_projection(m::GeoModel,mesh,dimension::Int,tag::Int,caller)
    input=_model_mesh_projection_input(m,mesh,dimension,tag,caller)
    input===nothing && return nothing
    authority=input.authority;actual=authority.mesh;stored=actual.entity_data
    closure=_model_mesh_projection_closure(m,dimension,tag,caller)
    included=falses(nnodes(actual));included[input.rows].=true
    for owner in closure
        product=get(authority.products,owner,nothing)
        product===nothing || (included[product.rows].=true)
    end
    rows=copy(input.rows)
    missing=Int32[row for row in eachindex(included) if included[row]]
    sort!(missing;by=row->authority.private_rank[row])
    present=falses(nnodes(actual));present[rows].=true
    for row in missing
        present[row] || push!(rows,row)
    end
    mapping=zeros(Int32,nnodes(actual));mapping[rows]=Int32.(1:length(rows))
    all(owner->owner in closure,stored.node_entities[rows]) || throw(ArgumentError(
        "$caller: actual projection references a node outside its entity closure"))
    blocks=ElementBlock[];owners=Vector{Int32}[];labels=Vector{UInt64}[]
    for (block_index,block) in enumerate(actual.blocks)
        selected=Int[column for column in axes(block.nodes,2) if
            msh_dimension(block.msh)<dimension &&
            (msh_dimension(block.msh),stored.block_entities[block_index][column]) in closure]
        isempty(selected) && continue
        push!(blocks,ElementBlock(block.msh,mapping[block.nodes[:,selected]],block.tags[selected]))
        push!(owners,stored.block_entities[block_index][selected])
        push!(labels,stored.external_element_tags[block_index][selected])
    end
    for (block_index,block) in enumerate(mesh.blocks)
        push!(blocks,ElementBlock(block.msh,block.nodes,block.tags))
        push!(owners,mesh.entity_data.block_entities[block_index])
        push!(labels,mesh.entity_data.external_element_tags[block_index])
    end
    entities=Dict{Tuple{Int,Int},MixedEntity}()
    names=Dict{Tuple{Int,Int},String}()
    for owner in closure
        entity=get(stored.entities,owner,nothing)
        entity!==nothing || throw(ArgumentError(
            "$caller: actual projection entity closure lacks current metadata"))
        expected=_model_projection_physical_tags(m,owner[1],Int(owner[2]))
        entity.physical_tags==expected || throw(ArgumentError(
            "$caller: actual projection physical memberships no longer match the model"))
        entities[owner]=entity
        for physical in expected
            key=(owner[1],Int(physical))
            haskey(m.physical_names,key) && (names[key]=m.physical_names[key])
        end
    end
    metadata=Elements._native_mixed_entity_data(entities;node_entities=stored.node_entities[rows],
        node_parametric=stored.node_parametric[rows],external_node_tags=stored.external_node_tags[rows],
        block_entities=owners,external_element_tags=labels)
    links=MixedPeriodicLink[]
    for link in actual.periodic_links
        (link.dim,link.slave_entity) in closure || continue
        (link.dim,link.master_entity) in closure || throw(ArgumentError(
            "$caller: actual periodic relation crosses the projection entity closure"))
        slave=mapping[link.slave_nodes];master=mapping[link.master_nodes]
        all(>(0),slave) && all(>(0),master) || throw(ArgumentError(
            "$caller: actual periodic relation references an absent projected node"))
        push!(links,MixedPeriodicLink(link.dim,link.slave_entity,link.master_entity,
            slave,master;affine=link.affine))
    end
    output=MixedMesh(actual.coords[:,rows],blocks;physical_names=names,entity_data=metadata,
        elementary_entities=owners,periodic_links=links)
    diagnostic=validate(output;reject_duplicate_point_cells=false)
    diagnostic.ok || throw(ArgumentError(
        "$caller: invalid current actual projection — "*join(diagnostic.messages,"; ")))
    return output
end
